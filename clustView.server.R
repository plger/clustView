#' clustView.server
#'
#' Server-side function of the `clustView` shiny app.
#'
#' @param seurat An object of class `seurat`, or a path to a RDS of such an object.
#' @param enableClustree Logical; whether to draw the clustering tree 
#' (default TRUE). Disable to increase speed.
#'
#' @export
clustView.server <- function( seurat,
                              enableClustree = T
                             ){
    library(shinycssloaders)
    library(data.table)
    library(DT)
    library(ggplot2)
    library(cowplot)
    
    function(input, output, session) {

        if(is.character(seurat)) seurat <- readRDS(seurat)
        
        # we make a local, reactive copy of the meta.data 
        # and cluster annotation, in order to be able to
        # update it
        sera <- reactiveValues( metadata=seurat@meta.data, 
                                anno=seurat@misc$clusterAnnotation )
        
        # preparing basic prefixes and resolutions
        cn <- colnames(seurat@meta.data)
        cn <- cn[grep("res\\.[0-9\\.]+$",cn)]
        allres <- t(sapply(cn,FUN=function(x){ strsplit(x,"res.",fixed=T)[[1]] }))
        allres <- split(allres[,2],paste0(allres[,1],"res."))
        prefs <- names(allres)
        
        if(is.null(seurat@misc$clusterAnnotation)) stop("Seurat object has no cluster annotation (in `seurat@misc$clusterAnnotation`)")
        
        prefs <- intersect(prefs, names(seurat@misc$clusterAnnotation))
        if(length(prefs)==0) stop("Empty cluster annotation, or mismatch with resolution prefixes from `seurat@meta.data`")
        allres <- allres[prefs]
        
        gocats <- unique(unlist(lapply(seurat@misc$clusterAnnotation, FUN=function(x){ names(x$go) })))
        
        updateSelectInput(session, 'prefix', choices = names(allres))
        updateSelectInput(session, 'resolution', choices = allres)
        updateSelectInput(session, 'space', choices = names(.getDimRed(seurat)))
        
        output$go_ui <- renderUI({
            lapply(gocats, w=max(c(4,floor(12/length(gocats)))), FUN=function(x,w){
                box( title=x, collapsible=T, width=w, 
                     withSpinner(DTOutput(paste0("go_",x)))
                   )
            })
        })

        observe({
            req(input$prefix)
            ch <- sort(as.character(allres[[input$prefix]]))
            updateSelectInput(session, 'resolution', choices=ch)
        })

        # use this to access the currently-selected identities
        idents <- reactive({
            req(input$resolution)
            sera$metadata[[paste0(input$prefix, input$resolution)]]
        })
        
        observeEvent(input$resolution, {
            ch <- sort(as.character(unique(idents())))
            updateSelectInput(session, 'cluster', choices = ch)
        })
        
        
        ###############
        ## Begin clustree tab
        
        ro_clustree <- reactive({
            if (!enableClustree)
                return(NULL)
            library(clustree)
            df <- cbind(rep("0", ncol(seurat)), sera$metadata)
            colnames(df)[1] <- paste0(input$prefix, "0")
            clustree(df, prefix = input$prefix, return = 'plot')
        })
        
        output$clustree_msg <- renderText({
            if (enableClustree)
                return(NULL)
            "Clustering tree disabled. Use `enableClustree=T` to enable."
        })
        
        output$clustree <- renderPlot({
            ro_clustree()
        })
        
        observeEvent(input$clustree_click, {
            p <- ro_clustree()
            if (is.null(p)) return(NULL)
            np <- nearPoints(p$data, input$clustree_click, maxpoints=1)
            if (!is.null(np) && length(np) > 0) {
                updateSelectInput(session, 'resolution', selected = paste0(input$prefix, np[[input$prefix]]))
                updateSelectInput(session, 'cluster', selected = idents()[np$cluster])
                updateTabItems(session, "tabs", selected = "details")
            }
        })
        
        
        ## End clustree tab
        ###############
        ## Begin overview tSNE
        
        output$tsne_overview <- renderPlot({
            Idents(seurat) <- idents()
            DimPlot(
                seurat,
                reduction = input$space,
                group.by = "ident",
                pt.size = 1.5,
                do.label = T,
                label.size = 7,
                do.return = T
            ) + theme_cowplot() 
        })
        
        observeEvent(input$overviewPlot_click, {
            w <- .getNearestPoint( input$overviewPlot_click,
                                   .getDimRed(seurat)[[input$space]]@cell.embeddings )
            updateSelectInput(session, 'cluster', selected = idents()[w])
            updateTabItems(session, "tabs", selected = "details")
        })
        
        ## End overview tSNE
        ###############
        ## Begin cluster details
        
        output$tsne_detail <- renderPlot({
            DimPlot(
                seurat,
                reduction = input$space,
                cells.highlight = list(selected = which(idents() == input$cluster)),
                sizes.highlight = 2,
                cols.highlight = "#3C8DBC",
                do.return = T,
                no.axes = T,
                no.legend = T
            ) + theme_cowplot() + guides(color = FALSE)
        })
        
        observeEvent(input$detailPlot_click, {
            w <- .getNearestPoint( input$detailPlot_click,
                                   .getDimRed(seurat)[[input$space]]@cell.embeddings )
            if (!is.null(w) && length(w) > 0){
                updateSelectInput(session, 'cluster', selected = idents()[w])
            }
        })
        
        for(x in gocats){
            output[[paste0("go_",x)]] <- renderDT({
                go <- sera$anno[[input$prefix]]$go[[x]][[input$resolution]][[input$cluster]]
                if(is.null(go)) return(NULL)
                datatable(go)
            })
        }

        output$markers <- renderTable({
            mrks <- sera$anno[[input$prefix]]$markers[[input$resolution]]
            mrks <- mrks[which(mrks$cluster == input$cluster),,drop=F]
            if(nrow(mrks) == 0) return(NULL)
            tmp <- t(sapply( mrks$gene, FUN = function(x) {
                strsplit(as.character(x), ".", fixed = T)[[1]]
            }))
            mrks$ensembl <- tmp[,1]
            mrks$symbol <- tmp[,2]
            mrks$pval <- format(as.numeric(mrks$p_val), digits=2, scientific=T)
            mrks[,c("ensembl","symbol","pval")]
        })
        
        output$rename_msg <- renderText({
            if(input$newname %in% levels(idents())){
                return( paste0(
                        "There is already a cluster named '",
                        input$newname,
                        "' in this clustering/resolution; ",
                        "no action will be taken." ) )
            }
            return(NULL)
        })
        observeEvent(input$save_newname, {
            newname <- input$newname
            lvls <- levels(idents())
            if( gsub(" ","",newname) != "" &&
                !(newname %in% lvls) ){
                var <- paste0(input$prefix, input$resolution)
                lvls[which(lvls==input$cluster)] <- newname
                levels(sera$metadata[[var]]) <<- lvls
                nn <- names(sera$anno[[input$prefix]]$markers)
                nn[which(nn==input$cluster)] <- newname
                names(sera$anno[[input$prefix]]$markers) <<- nn
                for(x in gocats){
                    nn <- names(sera$anno[[input$prefix]]$go[[gocats]][[input$resolution]])
                    nn[which(nn==input$cluster)] <- newname
                    names(sera$anno[[input$prefix]]$go[[gocats]][[input$resolution]]) <<- nn
                }
                updateTextInput(session, 'newname', value='')
                updateSelectInput(session, 'cluster', choices=lvls, selected=newname)
            }
        })
        
        ## End cluster details
        ###############
        ## Begin download
        
        output$downloadRDS <- downloadHandler(
            filename="seurat.rds",
            content = function(file) {
                seurat@misc$clusterAnnotation <- sera$anno
                seurat@meta.data <- sera$metadata
                Idents(seurat) <- idents()
                saveRDS(seurat, file)
            }
        )
        
    }
}

# returns the list of available DimRed; works with seurat v2 / v3
.getDimRed <- function(seurat){
    if ("dr" %in% slotNames(seurat)) return(seurat@dr)
    return(seurat@reductions)
}

.getNearestPoint <- function(event, a){
    which.min(colSums(abs(t(a[,1:2])-c(event$x,event$y))))
}
