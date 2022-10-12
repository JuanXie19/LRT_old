Il#' cluster inferred clonotype trajectories
#' @param Combined a \code{CombinedDataSet} object
#' @param trajectory a \code{TrajectoryDataSet} object
#'
#' @export
#' @import dplyr
#' @import shiny
#' @import shinyBS
#' @importFrom dtwclust tsclust
#' @import ggplot2
#' @import magrittr
#' @import ggpubr
#' @import cluster
#' @import parallel
#' @import princurve
#' @import tradeSeq
#' @import Seurat
#' @import pheatmap
#' @import shinysky
#' @import sbm
#' @import SingleCellExperiment
#' @examples
#' TCR <-read.csv("/PATH/TO/YOUR/scTCR-seqData/",header=T)
#' load('Mice.sub2.rda')
#' Combined <- getCombinedDataSet(TCR,Mice.sub)
#' Trajectory <- getClonotypeLineages(Combined,start.clus = NULL, end.clus = NULL, dist.method = 'simple', use.median = TRUE)
#' shinyClust(Trajectory,Combined)


shinyClust <- function(TrajectoryDataSet = NULL, CombinedDataSet = NULL){


library(SingleCellExperiment)
ui <- fluidPage(

  # App title ----
  titlePanel("shinyClust"),

  fluidRow(
	 column(2,
		wellPanel(	
			actionButton('btn_go','Run clustering'),
			hr(),
			bsCollapse(id="clusterInfo", open='Clustering options',

                   bsCollapsePanel("Clustering options",
								radioButtons(inputId = "clustertype",
                                               label = "Clustering algorithm",
                                               choices = c("SBM","Partitional", "Hierarchical"),
                                               selected = 'SBM'),
								uiOutput("clusterdetails")							
								)                                  
                   )					   			
				)),
		column(10,
			
		tabsetPanel(
			tabPanel('Clustering exploration',
			fluidRow(
				column(6,
					uiOutput('sbmplot'),
					uiOutput('sbmtext')
				),
				column(6,
					uiOutput('silhouttePlot'),
					uiOutput('silhoutte_text')
				))),
			tabPanel('Trajectory clustering',
			fluidRow(
				selectInput(inputId = 'clusterIndex',
							label = 'Choose a cluster to display',
							choices = c('cluster 1','cluster 2' ),
							selected = 'cluster 1')
				),
			fluidRow(
				column(6,
				h4('Trajectories'),
				uiOutput(outputId = 'plotLin'),
				uiOutput('downloadLin'),
				helpText('This plot shows the trajectories in a selected cluster.',
					'Back lines represent all trajectories.', 'Dots represent all cells')),

				column(6,
				h4('Cell distribution'),
				uiOutput(outputId = 'plotUMAP'),
				uiOutput('downloadUMAP'),
				helpText('This plot shows the distribution of cells on the UMAP.', 
				'Solid triganular points are cells invovled in this cluster of lineages,', 
				'circles represent cells not in this cluster' ))
		
		
			)
		),
		tabPanel('DE analysis',
			fluidRow(
				fluidRow(
					column(3,selectInput(inputId = 'clusterIndex1',
							label = 'Choose a cluster to predict DE genes',
							choices = c('cluster 1','cluster 2' ),
							selected = 'cluster 1')),
					column(3, radioButtons(inputId='filterGene',
									  label = 'Filter gene by',
									  choices = c('adjusted p value','log fold change' ),									
									  selected = 'adjusted p value',inline=TRUE)),
					column(3,uiOutput("filterdetails")),
					column(3,radioButtons(inputId = 'dePlot',
										label = 'Visualization',
										choices = c('model fit','heatmap'),
										selected = 'heatmap',inline=TRUE)
							)),
				
				fluidRow(
					column(6,
					busyIndicator(text = "Predicting DE genes, which may take 5~10 minutes, please wait..."),
					DT::dataTableOutput('deGenetable')
					),
					column(6,plotOutput('plotDE'))
					))
			),
	
		))))	   



	server <- function(input,output,session) {
	
		output$clusterdetails <- renderUI({
			if(is.null(input$clustertype))
				return()
				
			switch(input$clustertype,
				'SBM' = list(radioButtons(inputId = 'sbmCluster',
										label = 'Way to decide number of clusters',
										choices = c('data driven','user specified'),
										selected = 'data driven'),
							uiOutput('sbmdecision')
							),
				'Hierarchical' = list(
									numericInput(inputId = 'clustNumH',
												label = 'Desired number of clusters',
												min = 2, value =6),
									numericInput(inputId = 'setSeedH',
											label = 'set seed',value=123456),
									selectInput(inputId = 'lineage',
												label = 'hc.lineage',
												choices = c('complete','average','single','ward.D'),
												selected = 'average')
									 ),
				'Partitional' = list(
									numericInput(inputId = 'clustNumP',
												label = 'Desired number of clusters',
												min = 2, value =6),
									numericInput(inputId = 'setSeedP',
											label = 'set seed',value=123456),
									numericInput(inputId = 'iterMax',
												label = 'iterMax',
												min = 10, max=150, value = 100),
										selectInput(inputId = 'centroid',
											    label = 'centroid',
												choices = c('pam','dba','mean','median'),
												selected = 'pam')))			
		})
		output$sbmdecision <- renderUI({
			switch(input$sbmCluster,
				'user specified' = numericInput(inputId = 'clustNumS',
												label = 'Desired number of clusters',
												min = 2, value =6)
			)
		})

		df <- TrajectoryDataSet@lineages
		UMAP1 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,1]
		UMAP2 <- CombinedDataSet@seurat@reductions$umap@cell.embeddings[,2]
		cdr3 <- CombinedDataSet@seurat$cdr3
		cell.type <- CombinedDataSet@seurat$clusters
		G.df <- data.frame(UMAP1,UMAP2,cell.type,cdr3)

		itermax <- reactive(input$iterMax)
		centroid <- reactive(input$centroid)
		seedH <- reactive(input$setSeedH)
		seedP <- reactive(input$setSeedP)
		kP <- reactive(input$clustNumP)
		kH <- reactive(input$clustNumH)
		kS <- reactive(input$clustNumS)
		

		method <- reactive(input$lineage)
		#CHOICE <-reactive(input$radio)
		
		#### run SBM as default
		d <- dtwclust::tsclust(df,type='h',k=6,seed=123456,control=hierarchical_control(method='average'))@distmat  # DTW distance matrix
		## convert distance matrix to adjacency matrix based on KNN
		adjact <- matrix(0,nrow=nrow(d),ncol=nrow(d))
		colnames(adjact) <-rownames(adjact) <-rownames(d)
		diag(adjact) <- 1
		n <- ceiling(sqrt(length(df)))
		for (i in 1:nrow(d)){
			nn <- names(sort(d[i,])[1:n])
			INDEX <-which(colnames(adjact)%in% nn)
			adjact[i,INDEX] <- 1
		}
		
		
## run lineage clustering
	
		mySBM <- eventReactive(input$btn_go,{
				adjact %>% 
				estimateSimpleSBM("bernoulli", dimLabels = 'lineage', estimOptions = list(verbosity = 0, plot = FALSE))		
			})
		
		

		observeEvent(mySBM(),{
			updateSelectInput(session,'clusterIndex',choices = paste('cluster',1:mySBM()$nbBlocks))
		})
		
		observeEvent(mySBM(),{
			updateSelectInput(session,'clusterIndex1',choices = paste('cluster',1:mySBM()$nbBlocks))
		})

		
		
		clusterRST <- eventReactive(input$btn_go,{
				if (input$sbmCluster =='data driven')
					x <- switch(input$clustertype,
						'SBM' = mySBM()$memberships,
						'Hierarchical' = {
									temp <- tsclust(df,type='h',k=kH(),seed=seedH(),control=hierarchical_control(method=method()))
									temp@cluster},
						'Partitional' = {
									temp <- tsclust(df,type='p',k=kP(),seed=seedP(),iter.max=itermax(),centroid=centroid())
									temp@cluster})	
				if (input$sbmCluster == 'user specified')
					x <- switch(input$clustertype,
						'SBM' = {mySBM()$setModel(kS())
								 mySBM()$memberships
								},
					'Hierarchical' = {
									temp <- tsclust(df,type='h',k=kH(),seed=seedH(),control=hierarchical_control(method=method()))
									temp@cluster},
						'Partitional' = {
									temp <- tsclust(df,type='p',k=kP(),seed=seedP(),iter.max=itermax(),centroid=centroid())
									temp@cluster})
				return(x)
						
		})
		

		observeEvent(clusterRST(),{
			updateSelectInput(session,'clusterIndex',choices = paste('cluster',1:max(clusterRST())))
		})
		
		observeEvent(clusterRST(),{
			updateSelectInput(session,'clusterIndex1',choices = paste('cluster',1:max(clusterRST())))
		})
	
		
		clusterID <- reactive(as.numeric(strsplit(input$clusterIndex, ' ')[[1]][2]))  # used on clustering analysis
		clusterID1 <- reactive(as.numeric(strsplit(input$clusterIndex1, ' ')[[1]][2]))  # used on DE analysis
					
		UMAP <- function(CLONES){		
			#CLONES <- names(df)[which(mySBM()$memberships==clusterID())]
			INDEX <- which(G.df$cdr3 %in% CLONES)
			G.df.highlight <- G.df[INDEX,]
			ggplot(G.df,aes(UMAP1,UMAP2,color=cell.type))+geom_point(alpha=0.3,shape=1,size=0.5)+
				geom_point(data=G.df.highlight,aes(UMAP1,UMAP2),size=1.2,shape=17)+theme_bw()		
		}
		
		observeEvent(input$btn_go,{
			output$plotUMAP <- renderUI({
				plotOutput('sbmumap')
			})
			
		})
		
		
		
		output$sbmumap <- renderPlot({
			CLONES <- names(df)[which(clusterRST()==clusterID())]
			print(UMAP(CLONES))
		})
		

		
		observeEvent(input$btn_go,{
			output$downloadUMAP <-renderUI({
				downloadButton('downloadsbmumap',label = 'Download plot')
			})
		})
		

		output$downloadsbmumap <- downloadHandler(
			filename = function(){
				paste(input$clusterIndex,'_plot2.pdf',sep='')},
				content = function(file){
					pdf(file)
					CLONES <- names(df)[which(clusterRST()==clusterID())]
					print(UMAP(CLONES))
					dev.off()
				}			
		)
		
		## plot for default sbm clustering
		LIN <- function(INDEX){
			
			lin <- df[INDEX]  # the lineages in a particular cluster
			DF   <- cbind(cdr3=rep(names(lin),sapply(lin,nrow)),do.call(rbind,lin))
			colnames(DF) <-c('cdr3','UMAP1','UMAP2')
			
			ggplot(NULL, aes(UMAP1,UMAP2)) +geom_point(data=G.df,aes(col=cell.type),size=0.2,shape=19,alpha=0.2)+ 
					geom_line(data=DF,aes(group=cdr3))+theme_bw()						
		}
		
		observeEvent(input$btn_go,{
			output$plotLin <- renderUI({
				plotOutput('sbmlin')
			})
			
		})
		
		
		output$sbmlin <- renderPlot({
			INDEX <- which(clusterRST() == clusterID())
			print(LIN(INDEX))
		})
		

		
		## download Lineage plots
		observeEvent(input$btn_go,{
			output$downloadLin <-renderUI({
				downloadButton('downloadsbmlin',label = 'Download plot')
			})
		})
		

		
		output$downloadsbmlin <- downloadHandler(
			filename = function(){
				paste(input$clusterIndex,'_plot2.pdf',sep='')},
				content = function(file){
					pdf(file)
					#INDEX <- which(mySBM()$memberships == clusterID())
					INDEX <- which(clusterRST() == clusterID())
					print(LIN(INDEX))
					dev.off()
				}			
		)
		
		
###### DE analysis
		output$filterdetails <- renderUI({
			if(is.null(input$filterGene))
				return()
				
			switch(input$filterGene,
				"adjusted p value" = numericInput(inputId = 'pvalueThresh',
												label = 'p.adjust threshold',
												min = 0.0001, max =0.5,value = 0.05),
												
				"log fold change" = numericInput(inputId = 'logFC',
												label = 'logFC threshold',
												min = 0.1, value = 1))												
		})
		
		META <- CombinedDataSet@seurat@meta.data
		cell.umap <- CombinedDataSet@seurat@reductions$umap@cell.embeddings
		DATA <- CombinedDataSet@seurat@assays$RNA  # read counts
		
		
		stored <- reactiveValues(CLONES=NULL)
		
		if(FALSE){
		observeEvent(c(
					input$btn_go,
					input$clusterIndex1),
					{
					stored$CLONES <- gsub('[[:digit:]]+', '',names(df)[(which(mySBM()$memberships == clusterID1()))])
		})
		}
		observeEvent(c(
					input$btn_go,
					input$clusterIndex1),
					{
					stored$CLONES <- gsub('[[:digit:]]+', '',names(df)[(which(clusterRST() == clusterID1()))])

		})
		
		
		gamFIT <- reactive({
				CLONES <- stored$CLONES ## clones that belong to the lineage cluster
				Cluster.df <- META[which(META$cdr3%in%CLONES),]  # the meta data for those clones
				cell.umap.sub <- cell.umap[which(rownames(cell.umap)%in% rownames(Cluster.df)),]  # the UMAP for them
				fit <- principal_curve(cell.umap.sub,smoother='lowess')
				DATA.sub <- DATA[,which(colnames(DATA)%in% rownames(cell.umap.sub))]  # read counts for cells belonging to the lineage cluster
				DATA.sub <-DATA.sub[which(rownames(DATA.sub) %in% VariableFeatures(CombinedDataSet@seurat)),]  # just use highly variable genes
				weights <- as.matrix(rep(1,dim(DATA.sub)[2]))
				pseudotime <- as.matrix(fit$lambda)
				sce <- fitGAM(DATA.sub,pseudotime = pseudotime,cellWeights = weights,nknots=6,verbose=FALSE,parallel = TRUE)   # takes a long time
				colData(sce)$clusters <- Cluster.df$clusters
				return(sce)
		})
		
		assoRes <- reactive(associationTest(gamFIT()))
		#colnames(assoRes()) <- c('waldStat','raw p value','meanLogFC')
		#FDR <- reactive(round(p.adjust(assoRes()$pvalue,method="bonferroni"),3))
		FDR <- reactive(p.adjust(assoRes()$pvalue,method="bonferroni"))

		assoRes1 <- reactive({
			data.frame(round(assoRes(),3),adjusted.pvalue=round(FDR(),3))
		})
		
		filtered.sce <- reactive({ 
			switch(input$filterGene,
					'adjusted p value' = assoRes1() %>% filter (adjusted.pvalue <= input$pvalueThresh),
				    'log fold change' = assoRes1()%>% filter(meanLogFC >= input$logFC))
		})
		
		## may try waiter https://waiter.john-coene.com/#/hostess
		output$deGenetable = DT::renderDataTable(
							filtered.sce(),
							 server=FALSE,extensions = 'Buttons',options = list(
										order = list(5,'asc'),dom = 'frtipB',buttons = list('copy','print',list(extend = 'collection',buttons = c('csv','excel'),text='Download table'))))
										
		
		s <- reactive({
			req(input$deGenetable_rows_selected)
			input$deGenetable_rows_selected
			})
			
		smoothPlot <- function(){		
			CLONES <- stored$CLONES # clones that belong to the lineage cluster
			Cluster.df <- META[which(META$cdr3%in%CLONES),]  # the meta data for those clones
			cell.umap.sub <- cell.umap[which(rownames(cell.umap)%in% rownames(Cluster.df)),]  # the UMAP for them
			DATA.sub <- DATA[,which(colnames(DATA)%in% rownames(cell.umap.sub))]  # read counts for cells belonging to the lineage cluster
			DATA.sub <-DATA.sub[which(rownames(DATA.sub) %in% VariableFeatures(CombinedDataSet@seurat)),]  # just use highly variable genes
			#remove(Cluster.df)
			
			#oStart <- order(filtered.sce()$waldStat, decreasing = TRUE)
			#sigGeneStart <- names(gamFIT())[oStart]
			plotSmoothers(gamFIT(), DATA.sub, gene = rownames(filtered.sce())[s()],pointCol='clusters')
			
		}
		
		seurat <- CombinedDataSet@seurat
		## heatmap
		HEATMAP <- function(){
				CLONES <- stored$CLONES
				Cluster.df <- META[which(META$cdr3%in%CLONES),]
				seurat.filtered <- seurat[,colnames(seurat)%in%Cluster.df$barcode]
				matrix_to_plot <- AverageExpression(seurat.filtered,assay='SCT',slot='data',group_by='clusters')
				matrix_to_plot <- matrix_to_plot[[1]]
				matrix_to_plot <- matrix_to_plot[rownames(filtered.sce()), ]

				col_df <- data.frame(row.names = colnames(matrix_to_plot), cluster = colnames(matrix_to_plot))
				pmap <-pheatmap(matrix_to_plot, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE, show_colnames = TRUE, show_rownames = FALSE,
					annotation_col = col_df)
				print(pmap)
		}
		
		output$plotDE <- renderPlot({
			switch(input$dePlot,
				"model fit" = print(smoothPlot()),
				"heatmap" = print(HEATMAP())
	
			)

		})
				
	### clustering exploration functions	

		observeEvent(input$btn_go,{
			output$sbmplot <- renderUI({
				switch(input$clustertype,
					'SBM' = plotOutput('sbmICL'),
					'Hierarchical' = plotOutput('MDS'),
					'Partitional' = plotOutput('MDS')
				)
			})		
		})
		
		output$sbmICL <- renderPlot({
			mySBM()$storedModels %>%  ggplot() + aes(x = nbBlocks, y = ICL)  + geom_line() + geom_point(alpha = 0.5)
		})
		
		observeEvent(input$btn_go,{
			output$silhouttePlot <- renderUI({
				if(input$sbmCluster == 'data driven')
				 x <- switch(input$clustertype,
					'SBM' = plotOutput('adjacencyplot'),
					'Hierarchical' = plotOutput('silhouette'),
					'Partitional' = plotOutput('silhouette'))
				if(input$sbmCluster == 'user specified')
					x <- switch(input$clustertype,
					'SBM' = plotOutput('adjacencyplot1'),
					'Hierarchical' = plotOutput('silhouette'),
					'Partitional' = plotOutput('silhouette'))
				return(x)
			})		
		})
		
		output$adjacencyplot <- renderPlot({
			plot(mySBM(), type = "data", dimLabels  = list(row = 'lineage', col= 'lineage'))

		})
		output$adjacencyplot1 <- renderPlot({
			mySBM()$setModel(kS())
			mySBM()$memberships
			plot(mySBM(), type = "data", dimLabels  = list(row = 'lineage', col= 'lineage'))

		})
		
		
		MDS <- function(){
			d <- dtwclust::tsclust(df,type='h',k=6,seed=123456,control=hierarchical_control(method='average'))  # DTW distance matrix
			DIST <- d@distmat
			t <- colnames(DIST)
			# Compute MDS
			mds <- DIST %>%
				cmdscale() %>%
					as_tibble()
			colnames(mds) <- c("Dim.1", "Dim.2")

			mds <- mds %>%
				mutate(groups = as.factor(clusterRST()))
			# Plot MDS
			ggscatter(mds, x = "Dim.1", y = "Dim.2",
				label = NULL,color = "groups",
				#palette = "jco",
				size = 1,ggtheme = theme_bw(),
				ellipse = TRUE,
				ellipse.type = "convex",
				repel = TRUE,ggrepel.max.overlaps = Inf)
		}
		output$MDS <- renderPlot({
			print(MDS())
		})
		
		Silhouette <- function(){
			d <- dtwclust::tsclust(df,type='h',k=6,seed=123456,control=hierarchical_control(method='average')) # DTW distance matrix
			sil <- cluster::silhouette(clusterRST(),dist= d@distmat,title=title(main = 'Good'))
			plot(sil,main='Silhouette plot')
		}
		
		output$silhouette <- renderPlot({
			print(Silhouette())
		})
		
		## helptext
	
		output$sbmICLtext <- renderText({
			paste("The plot shows how ICL changes with number of blocks(clusters). We'd prefer a model with higher ICL ")
		})
		output$sbmadjatext <- renderText({
			paste("The plot shows the adjacency matrix for selected model, which is reordered according to the memberships")

		})
		
		observeEvent(input$btn_go,{
			output$sbmtext <- renderUI({
				switch(input$clustertype,
					'SBM' = textOutput('sbmICLtext'),
					'Hierarchical' = textOutput('MDStext'),
					'Partitional' = textOutput('MDStext')
				)
			})		
		})
		
		observeEvent(input$btn_go,{
			output$silhoutte_text <- renderUI({
				switch(input$clustertype,
					'SBM' = textOutput('sbmadjatext'),
					'Hierarchical' = textOutput('sil_text'),
					'Partitional' = textOutput('sil_text')
				)
			})		
		})
		
		output$MDStext <- renderText({paste("This MDS plot visualize the distance between lineages in a two dimensional space")})
		output$sil_text <- renderText({paste("This plot shows silhouette scores(ranges from -1 to 1). High value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters")})
		

		
	}

	shinyApp(ui = ui, server = server)
}