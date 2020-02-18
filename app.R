#library(shiny)
library(ggplot2)
library(Seurat)

#Sweave('melton_Resources_VeresData_Stage5_02JAN20.Rnw')
load('cpobj5.rda')
load('cpobj6.rda')
load('usecolors.rda')


## UI
ui <- fluidPage(
    titlePanel(h3('Human Beta-Cell Single Cell Assay - Veres A, etal.)'), windowTitle='betqcelltypes'),
    sidebarLayout(
        sidebarPanel(
            conditionalPanel(
                condition = 'input.thetabs == 1',
                textInput(inputId='symbol', 'Symbol of Gene to map on clusters', value='INS'),
                verbatimTextOutput('info'))
        ),
        mainPanel(
            tabsetPanel(id='thetabs', selected=1, type='tabs',
                        tabPanel('Stage5', value=1, fluidRow(column(12, offset=1, plotOutput('mapping5', width='500px', height='500px')),
                                                            column(12, offset=1, plotOutput('stages5', width='500px', height='500px')))),
                        tabPanel('Stage6', value=1, fluidRow(column(12, offset=1, plotOutput('mapping6', width='500px', height='500px')),
                                                             column(12, offset=1, plotOutput('stages6', width='500px', height='500px')))),
                        tabPanel('Violin Plot', value=1, fluidRow(column(12, offset=1, plotOutput('vplot5', width='500px', height='400px')),
                                                             column(12, offset=1, plotOutput('vplot6', width='500px', height='400px'))))
                        
                        )
        )
    )
)

## Server
server <- function(input, output) {
    output$mapping5 <- renderPlot({
        output$info  <- renderText({'Only One Symbol'})

        symreact <- reactive({
            #validate(
             #   need(any(grepl(toupper(sub('\\.', '\\\\.', input$symbol)), ls(sym2name))) & input$symbol != '',
              #       "This symbol is not recognized or the gene is not expressed in any of the cells! Please, choose another gene.")
            #)
            toupper(input$symbol)
        })
    
        output$info <- renderText({symreact()})
        p <- FeaturePlot(cpobj5, features=symreact())
        ##p <- FeaturePlot(cpobj5, features='INS')
        plot(p)
    })

    output$stages5 <- renderPlot({
        straightcols <- unique(sub('[0-9]+', '', colors()))
        straightnogray <- straightcols[grep('grey|gray', straightcols, invert=TRUE)]
        #usecolors <- sample(straightnogray, 12)

        p <- DimPlot(cpobj5, group.by='Assigned_cluster', cols=usecolors, sizes.highlight=0.1, label=FALSE)
        p <- p +  ggtitle('Hs Islets Stage 5 Veres tSNE Coords') + theme(title=element_text(size=10))# + NoLegend()                                                     
        plot(p)
    })

    output$mapping6 <- renderPlot({
        output$info  <- renderText({'Only One Symbol'})

        symreact <- reactive({
            #validate(
             #   need(any(grepl(toupper(sub('\\.', '\\\\.', input$symbol)), ls(sym2name))) & input$symbol != '',
              #       "This symbol is not recognized or the gene is not expressed in any of the cells! Please, choose another gene.")
            #)
            toupper(input$symbol)
        })
    
        output$info <- renderText({symreact()})
        p <- FeaturePlot(cpobj6, features=symreact())
        plot(p)
    })

    output$stages6 <- renderPlot({
        straightcols <- unique(sub('[0-9]+', '', colors()))
        straightnogray <- straightcols[grep('grey|gray', straightcols, invert=TRUE)]
        #usecolors <- sample(straightnogray, 12)

        p <- DimPlot(cpobj6, group.by='Assigned_cluster', cols=usecolors, sizes.highlight=0.1, label=FALSE)
        p <- p +  ggtitle('Hs Islets Stage 6 Veres tSNE Coords') + theme(title=element_text(size=10))# + NoLegend()                                                     
        plot(p)
    })

    output$vplot5 <- renderPlot({
        output$info  <- renderText({'Only One Symbol'})

        symreact <- reactive({
            #validate(
             #   need(any(grepl(toupper(sub('\\.', '\\\\.', input$symbol)), ls(sym2name))) & input$symbol != '',
              #       "This symbol is not recognized or the gene is not expressed in any of the cells! Please, choose another gene.")
            #)
            toupper(input$symbol)
        })
    
        output$info <- renderText({symreact()})
        p <- VlnPlot(cpobj5, features=symreact(), slot = "counts", log = TRUE, group.by='Assigned_cluster', pt.size=0) + NoLegend()
        p <- p + NoLegend() + ggtitle(paste(symreact(), '(Stage5)'))
        plot(p)
     })

    output$vplot6 <- renderPlot({
        #output$info  <- renderText({'Only One Symbol'})

        symreact <- reactive({
            #validate(
             #   need(any(grepl(toupper(sub('\\.', '\\\\.', input$symbol)), ls(sym2name))) & input$symbol != '',
              #       "This symbol is not recognized or the gene is not expressed in any of the cells! Please, choose another gene.")
            #)
            toupper(input$symbol)
        })
    
        output$info <- renderText({symreact()})
        p <- VlnPlot(cpobj6, features=symreact(), slot = "counts", log = TRUE, group.by='Assigned_cluster', pt.size=0) + NoLegend()
        p <- p + NoLegend() + ggtitle(paste(symreact(), '(Stage6)'))
        plot(p)
    })

}

shinyApp(ui=ui, server=server)
    
