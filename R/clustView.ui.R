#' clustView.ui
#'
#' UI function of `clustView`.
#'
#' @export
clustView.ui <- function(){
    library(shiny)
    library(DT)
    library(shinycssloaders)
    library(shinydashboard)
    
    dashboardPage(
        dashboardHeader(title="clustViewer"),
        dashboardSidebar(
            sidebarMenu(id="tabs",
                menuItem("Clustree", tabName="tree"),
                menuItem("Clusters overview", tabName="overview"),
                menuItem("Cluster details", tabName="details"),
                menuItem("Download", tabName="download")
            ),
            selectInput('prefix', 'clustering', choices=c(), selectize=T),
            selectInput('resolution', 'resolution', choices=c(), selectize=T),
            selectInput('space', 'space', choices=c(), selectize=T)
        ),
        dashboardBody(
            tags$head(tags$style(type="text/css", '
                            .inlineDiv label { display: table-cell; vertical-align: middle; } 
                            .inlineDiv .form-group { display: table-row; }
                        ')),
            tabItems(
                tabItem("tree",
                        box(width=12, 
                            tags$div(style="font-weight: bold;", textOutput('clustree_msg')),
                            withSpinner(plotOutput("clustree", height='600px', click="clustree_click"))
                        )
                ),
                tabItem("overview",
                        box(width=12, 
                          withSpinner(plotOutput('tsne_overview', height='700px', click="overviewPlot_click"))
                        )
                ),
                tabItem("details",
                        fluidRow(
                            column( width=6, div(class="inlineDiv", selectInput('cluster', ' Cluster ', choices=c(), selectize=F) ) ),
                            column( width=6,
                             div(style="display: inline-block; vertical-align:top; width: 250px;", textInput('newname', label=NULL, placeholder = 'Enter new name')),
                             div(style="display: inline-block; vertical-align:top;", actionButton('save_newname', 'Rename cluster') ),
                             tags$p(textOutput('rename_msg'))
                           )
                        ),
                        box( withSpinner(plotOutput('tsne_detail', click="detailPlot_click")) ),
                        box( title = "markers", solidHeader=T, collapsible=T, 
                             div(style = 'height: 420px; overflow-x: scroll', tableOutput('markers') )
                        ),
                        box( title = "Prediction from dataset A" ), # not yet implemented
                        box( title = "Prediction from dataset B" ),
                        uiOutput("go_ui"),
                        div(style="clear: both;")
                ),
                tabItem("download",
                         box( tags$p("Download the Seurat object (with eventual modifications) in RDS format."),
                              downloadButton("downloadRDS", "Download RDS") )
                )
            ) # end tabItems
        )
    )
}