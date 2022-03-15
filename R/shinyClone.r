shinyClone <- function(CombinedDataSet =NULL){

# load libraries
	library(Seurat)
	library(dplyr)
# extract needed information
	UMAP1 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,1]
	UMAP2 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,2]
	cdr3 <- CombinedDataSet@seurat$cdr3
	G <- CombinedDataSet@seurat$clusters
	G.df <- data.frame(UMAP1,UMAP2,G,cdr3)
		
	CLONE <- unique(G.df$cdr3)
	n <-length(CLONE)
	
	ui <-fluidPage(
		titlePanel('Clonoal distribution on UMAP'),
		sidebarLayout(
		 sidebarPanel(
		 numericInput(inputId = 'cloneIndex',label = strong('Clonotype Id'),value = 2),
		 downloadButton('downloadPlot',label='Download Plot')
		 ),
		 mainPanel(
			textOutput(outputId = 'desc'),
			plotOutput(outputId = 'UMAP')
			
			
		 )
		)
	)
	
	server <- function(input,output){
		cloneInd <- reactive({
			req(input$cloneIndex)
			validate(need (input$cloneIndex>0 &input$cloneIndex <=length(CLONE),'Error:invalid range of clonotype index.The index should be greater than 0 and not larger than the number of unique clonotype'))
			validate(need(is.numeric(input$cloneIndex),'Error: cloneIndex should be a integer'))		
			INDEX <-which(G.df$cdr3 ==CLONE[input$cloneIndex])
			return(INDEX)
		})
		
		## description of the plotted clonotype
		output$desc <-renderText({
			G.df.highlight <- G.df[cloneInd(),]
			paste('This clonotype consists of ',length(cloneInd()),'cells','spanning',length(unique(G.df.highlight$G)),'cluster/clusters')
		})
		
		## plot the distribution
		plotUMAP <-function(){
			if(length(cloneInd())<3){
				print('too few cells')
			}
			
			G.df.highlight <- G.df[cloneInd(),]
			
			ggplot(G.df,aes(UMAP1,UMAP2,color=G))+geom_point(alpha=0.4,shape=1,size=0.6)+
				geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=2,shape=17)+
				theme_bw()+labs(title=CLONE[cloneInd()])
		
		}
		
		output$UMAP <-renderPlot({
			print(plotUMAP())
		})
		output$downloadPlot <-downloadHandler(
			filename = function(){
			paste(CLONE[cloneInd()],'_UMAP.pdf',sep='')},
			content = function(file){
			 pdf(file)
			 print(plotUMAP())
			 dev.off()
			}
		)
		
	}
	## create shiny object
	shinyApp(ui = ui, server = server)
}