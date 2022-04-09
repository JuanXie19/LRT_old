shinyClust <- function(TrajectoryDataSet = NULL, G.df = NULL){

# load libraries
	library(dplyr)
	library(shinyBS)
	library(shiny)
	library(dtwclust)

ui <- fluidPage(

  # App title ----
  titlePanel("Clustering analysis of clonotype lineages"),

  fluidRow(
	 column(2,
		wellPanel(	
			numericInput(inputId = 'clustNum',
						label = 'Desired number of clusters',
						min = 2, value =6),
			hr(),
			selectInput(inputId = 'clusterIndex',
							label = 'Choose a cluster to display',
							choices = c('cluster 1','cluster 2' ),
							selected = 'cluster 1'),
			hr(),
			radioButtons(inputId='radio',
									  label = 'Showing option for lineages in a cluster',
									  choices = list('show all'="1",'show representative'="2" ),									
									  selected = '1'),
			hr(),
			bsCollapse(id="clusterInfo", open='Clustering details',

                   bsCollapsePanel("Clustering details",
                                radioButtons(inputId = "clustertype",
                                               label = "Type of clustering",
                                               choices = c("Partitional", "Hierarchical"),
                                               selected = 'Hierarchical'),
									uiOutput("clusterdetails"),
									numericInput(inputId = 'setSeed',
												label = 'set seed',value=123456)
								)                                  
                   )
				   			
				)),
		column(5,
				h4('Lineages'),
				plotOutput('plotLin'),
				downloadButton('downloadPlot1',label = 'Download plot1'),
				helpText('This plot shows the lineages in a selected cluster.',
					'Back lines represent all lineages,', 'red one denotes the reprentative lineage. Dots represent all cells')),

		column(5,
				h4('Cell distribution'),
				plotOutput('plotUMAP'),
				downloadButton('downloadPlot2',label = 'Download plot2'),
				helpText('This plot shows the distribution of cells on the UMAP.', 
				'Solid triganular points are cells invovled in this cluster of lineages,', 
				'circles represent cells not in this cluster' ))

	)
)
		   



	server <- function(input,output,session) {
	
		output$clusterdetails <- renderUI({
			if(is.null(input$clustertype))
				return()
				
			switch(input$clustertype,
				'Hierarchical' = selectInput(inputId = 'lineage',
												label = 'hc.lineage',
												choices = c('complete','average','single','ward.D'),
												selected = 'average'),
				'Partitional' = list(numericInput(inputId = 'iterMax',
												label = 'iterMax',
												min = 5, max=100, value = 20),
										selectInput(inputId = 'centroid',
											    label = 'centroid',
												choices = c('pam','dba','mean','median'),
												selected = 'pam')))			
		})

		df <- TrajectoryDataSet@lineages

		itermax <- reactive(input$iterMax)
		centroid <- reactive(input$centroid)
		seed <- reactive(input$setSeed)
		k <- reactive(input$clustNum)

		method <- reactive(input$lineage)
		CHOICE <-reactive(input$radio)
		
		
		observeEvent(input$clustNum,{
			updateSelectInput(session,'clusterIndex',choices = paste('cluster',1:input$clustNum))
		})
		
		clusterID <- reactive(as.numeric(strsplit(input$clusterIndex, ' ')[[1]][2]))



		UMAP <- function(){
		## clustering 
			clusterRST <- reactive({
				switch(input$clustertype,
					'Hierarchical' = tsclust(df,type='h',k=k(),seed=seed(),control=hierarchical_control(method=method())),
					'Partitional' = tsclust(df,type='p',k=k(),seed=seed(),iter.max=itermax(),centroid=centroid())
			
				)
			})
			
		CLONES <- names(df)[(which(clusterRST()@cluster==clusterID()))]
		INDEX <- which(G.df$cdr3 %in% CLONES)
		G.df.highlight <- G.df[INDEX,]
		ggplot(G.df,aes(UMAP1,UMAP2,color=cell.type))+geom_point(alpha=0.3,shape=1,size=0.5)+
				geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=1.2,shape=17)+theme_bw()		
		}
		
		output$plotUMAP <- renderPlot({
			print(UMAP())
		
		})
		
		output$downloadPlot2 <- downloadHandler(
			filename = function(){
			paste(input$clusterIndex,'_plot2.pdf',sep='')},
			content = function(file){
				pdf(file)
				print(UMAP())
				dev.off()
			}		
		)
		
		
		LIN <- function(){
			clusterRST <- reactive({
				switch(input$clustertype,
					'Hierarchical' = tsclust(df,type='h',k=k(),seed=seed(),control=hierarchical_control(method=method())),
					'Partitional' = tsclust(df,type='p',k=k(),seed=seed(),iter.max=itermax(),centroid=centroid())
			
				)
			})
			
			INDEX <- which(clusterRST()@cluster == clusterID())
			
			CENTROID <- as.data.frame(clusterRST()@centroids[[clusterID()]])
			colnames(CENTROID) <- c('UMAP1','UMAP2')
			
			lin <- df[INDEX]  # the lineages in a particular cluster
			DF   <- cbind(cdr3=rep(names(lin),sapply(lin,nrow)),do.call(rbind,lin))
			colnames(DF) <-c('cdr3','UMAP1','UMAP2')
			
			switch(input$radio,
				"1" = ggplot(NULL, aes(UMAP1,UMAP2)) +geom_point(data=G.df,aes(col=cell.type),size=0.2,shape=19,alpha=0.2)+ 
					geom_line(data=DF,aes(group=cdr3))+theme_bw(),
				"2" = ggplot(NULL, aes(UMAP1,UMAP2)) +geom_point(data=G.df,aes(col=cell.type),size=0.2,shape=19,alpha=0.2)+ 
					geom_line(data=CENTROID,col='red',size=2)+theme_bw()
			)
			

		}
		
		
		output$plotLin <- renderPlot({
			print(LIN())
		
		})
		
		output$downloadPlot1 <- downloadHandler(
			filename = function(){
				paste(input$clusterIndex,'_plot1.pdf',sep='')
			},
			content = function(file){
				pdf(file)
				print(LIN())
				dev.off()
			
			}
		)
		
	}

	shinyApp(ui = ui, server = server)
}