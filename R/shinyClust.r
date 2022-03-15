shinyClust <- function(ClusterInput =NULL){

# load libraries
	library(dplyr)
	library(shinyBS)
	library(shiny)
	library(dtwclust)

ui <- fluidPage(

  # App title ----
  titlePanel("Clustering analysis of clonotype lineages"),

  sidebarLayout(
	wellPanel(	
	
		bsCollapse(id="clusterInfo", open="Clustering Details",
                   
                   bsCollapsePanel("Clustering Details",
                                   
                                   radioButtons(inputId = "clustertype",
                                               label = "Type of clustering",
                                               choices = c("Partitional", "Hierarchical"),
                                               selected = 'Hierarchical'),
									sliderInput(inputId = 'clustNum',
												label = 'Number of clusters',
												min = 2,max = 8,value = 2),
							
									## only show this if clustering type is partitional
									conditionalPanel(
										condition = "input.clustertype=='Partitional'",
										numericInput(inputId = 'iterMax',
												label = 'iterMax',
												min = 5, max=100, value = 20),
										selectInput(inputId = 'centroid',
											    label = 'centroid',
												choices = c('pam','dba','mean','median'),
												selected = 'pam'),
										numericInput(inputId = 'setSeed',
												label = 'seed',
												value = 123456)
									),							
									
									## only show this if clustering type is hierarchicha
									conditionalPanel(
										condition = "input.clustertype=='Hierarchical'",							
										selectInput(inputId = 'lineage',
												label = 'hc.lineage',
												choices = c('complet','average','single','ward.D'),
												selected = 'average')
									)									
								)                                  
                   )
				),
	mainPanel(
		tabsetPanel(
			tabPanel('Plot all lineages',
				plotOutput('plotLineages'),
				downloadButton('downloadLineage',label = 'Download lineage plot')
			),
			tabPanel('Plot clusters',
				conditionalPanel(
				condition = "input.clustertype=='Hierarchical'",
				plotOutput('plot_hc')),
				conditionalPanel(
				condition = "input.clustertype=='Partitional'",
				plotOutput('plot_p'))
			),
			tabPanel('Summary',
			conditionalPanel(
			condition = "input.clustertype=='Hierarchical'",
			verbatimTextOutput('summaryHC')),
			conditionalPanel(
			condition="input.clustertype=='Partitional'",
			verbatimTextOutput('summaryP'))
			),
			tabPanel('Cluster label',
			conditionalPanel(
			condition = "input.clustertype=='Hierarchical'",
			verbatimTextOutput('labelHC'),
			downloadButton('downloadTableHC',label = 'Download lineage cluster label')),
			conditionalPanel(
			condition="input.clustertype=='Partitional'",
			verbatimTextOutput('labelP'),
			downloadButton('downloadTableP',label = 'Download lineage cluster labels'))		
			)
		)
		
    )		   
  )
 )


	server <- function(input,output) {
		
		df <-ClusterInput@lineages
		
		itermax <-reactive(input$itermax)
		centroid <-reactive(input$centroid)
		seed <-reactive(input$setSeed)
		k <-reactive(input$clustNum)
	
	
		method <-reactive(input$lineage)
        ## plot all lineages
		
		plotLineages <-function(){
			UMAP1_min <- floor(min(sapply(df,function(x) min(x$UMAP_1))))-1
			UMAP1_max <- ceiling(max(sapply(df,function(x) max(x$UMAP_2))))+1
			UMAP2_min <- floor(min(sapply(df,function(x) min(x$UMAP_2))))-1
			UMAP2_max <- ceiling(max(sapply(df,function(x) max(x$UMAP_2))))+1
		
			plot(df[[1]]$UMAP_1,df[[1]]$UMAP_2,type='l',xlab='UMAP1',col='blue',ylab='UMAP2',xlim=c(UMAP1_min,UMAP1_max),ylim=c(UMAP2_min,UMAP2_max),pch=19,lwd=1.5)
			text(df[[1]]$UMAP_1,df[[1]]$UMAP_2,labels=rownames(df[[1]]),cex=0.5)
	
			for (i in 2:length(df)){
				lines(df[[i]]$UMAP_1,df[[i]]$UMAP_2,type='l',col="blue",pch=19,cex=0.5,lwd=1.5)
				text(df[[i]]$UMAP_1,df[[i]]$UMAP_2,labels=rownames(df[[i]]),cex=0.5)
			}	
		
		
		}
		output$plotLineages <-renderPlot({
			print(plotLineages())
		})
		
		output$downloadLineage <-downloadHandler (
			filename = function(){
				'AllClonotype Lineages.pdf'
			},
			content = function(file){
				pdf(file)
				print(plotLineages())
				dev.off()			
			}
		)
		
		## plot clusters of lineages
		output$plot_p <- renderPlot({
		
		  part.dtw <-tsclust(df,type='p',k=k(),seed=12345,iter.max=itermax(),centroid=centroid())
		  plot(part.dtw)
		})
	
		output$plot_hc <-renderPlot({
			hc.dtw <-tsclust(df,type='h',k=k(),control=hierarchical_control(method=method()))
			plot(hc.dtw,type='series')
		})
	
		output$summaryP <-renderPrint({
	
			part.dtw <-tsclust(df,type='p',k=k(),seed=12345,iter.max=itermax(),centroid=centroid())
			part.dtw
		})
	
		output$summaryHC <-renderPrint({
			hc.dtw <-tsclust(df,type='h',k=k(),control=hierarchical_control(method=method()))
			hc.dtw
	
		})
	
		output$labelP <- renderPrint({
	
			part.dtw <- tsclust(df,type='p',k=k(),seed=12345,iter.max=itermax(),centroid=centroid())
			part.dtw@cluster
		})
	
		output$labelHC <- renderPrint({
			hc.dtw <- tsclust(df,type='h',k=k(),control=hierarchical_control(method=method()))
			hc.dtw@cluster
	
		})
		
		output$downloadTableHC <- downloadHandler(
			filename = function(){
				'lineage cluster labels.csv'
			},
			content = function(file){
			hc.dtw <- tsclust(df,type='h',k=k(),control=hierarchical_control(method=method()))
			temp <- as.data.frame(hc.dtw@cluster)
			write.csv(temp,file)			
			}		
		)
		
		output$downloadTableP <- downloadHandler(
			filename = function(){
				'lineage cluster labels.csv'
			},
			content = function(file){
			part.dtw <- tsclust(df,type='p',k=k(),seed=12345,iter.max=itermax(),centroid=centroid())
			temp <- as.data.frame(part.dtw@cluster)
			write.csv(temp,file)			
			}		
		)
		
		
	}
	
	shinyApp(ui = ui, server = server)
}