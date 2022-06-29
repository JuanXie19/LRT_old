#' display the distribution of clonotype on UMAP
#' @param df a dataframe obtained from \code{getShinyinput}
#' @return shiny app
#' @export
#' @import Seurat
#' @import dplyr
#' @import DT
#' @import ggplot2

shinyClone <- function(G.df =NULL){


# extract needed information
  if(FALSE){
	UMAP1 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,1]
	UMAP2 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,2]
	cdr3 <- CombinedDataSet@seurat$cdr3
	G <- CombinedDataSet@seurat$clusters
	G.df <- data.frame(UMAP1,UMAP2,G,cdr3)
	}
	
	CLONES <- unique(G.df$cdr3)
	n <-length(CLONES)
	
	## find each clonotype involve how many cells
	
	cdr3.df <- data.frame(rbind(table(G.df$cdr3,G.df$cell.type)))
	cdr3.df$cell.num <- rowSums(cdr3.df)
	cdr3.df$cluster.num <-apply(cdr3.df,1, function(x) sum(x>0)-1) # count how many clusters it spans
	
	library(shiny)



	ui <-fluidPage(
		fluidRow(
			column(7,
				DT::dataTableOutput('clonetable')),
			column(5,
				textOutput(outputId = 'desc'),
				plotOutput('UMAP'),
				downloadButton('downloadPlot',label = 'Download Plot'))
		)
	)

server <-function(input,output,session){

  output$clonetable = DT::renderDataTable(cdr3.df,server=FALSE,extensions = 'Buttons',options = list(
										dom = 'Bfrtip',buttons = list('copy','print',list(extend = 'collection',buttons = c('csv','excel'),text='Download table')))
  )
  
  
  observe({
    req(input$clonetable_rows_selected)
    s = input$clonetable_rows_selected
	CLONE <- rownames(cdr3.df)[s]
	INDEX <- which(G.df$cdr3==CLONE)
	G.df.highlight <- G.df[INDEX,]
	## description of the plotted clonotype
	output$desc <-renderText({	
		paste('This clonotype consists of ',length(INDEX),'cells', 'spanning',length(unique(G.df.highlight$cell.type)),'cluster/clusters')	
	})
	
	plotUMAP <- function(){
		  ggplot(G.df,aes(UMAP1,UMAP2,color=cell.type))+geom_point(alpha=0.4,shape=1,size=0.6)+
          geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=2,shape=17)+
          theme_bw()+labs(title=CLONE)	
	}
	
	
	
    output$UMAP <-renderPlot({
      if (length(s)){
		print(plotUMAP())
	  }  
    })
	
	output$downloadPlot <-downloadHandler(
		filename = function(){
		paste(CLONE,'_UMAP.pdf',sep='')},
		content = function(file){
			 pdf(file)
			 print(plotUMAP())
			 dev.off()	
	})
    
  })
  
  
  
}

## create shiny object
shinyApp(ui = ui, server = server)
}
