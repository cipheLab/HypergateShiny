library("flowCore")
library("hypergate")
library("sp")
library("shinyjs")

logiclTransformCiphe <- function(flow.frame)
{

	#######################################################################################################
	#
	#
	#
	#
	#
	#######################################################################################################

	# no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
	# markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
	markers.transform <- colnames(flow.frame@description[["SPILL"]])

	list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
	list.index <- gsub("N","", list.index)
	list.index <- gsub("\\$P","", list.index)
		
	if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]]))
	{
		r.values <- unlist(lapply(list.index, function(x) 
			as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
		)	
	} else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]]))
	{
		r.values <- unlist(lapply(list.index, function(x) 
			as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
		)	
	} else 
	{
		r.values <- rep(90, length(list.index))
	}
	
	w.values <- (4.5-log10(262143/abs(r.values)))/2
	w.values[which(w.values<0)] <- 0.5
	w.values[which(is.infinite(w.values))] <- 0.5

	for(t in 1:length(markers.transform)){
		lgcl <- logicleTransform(w=w.values[t])
		flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
	}

	return(flow.frame)
}



server <- function(input, output, session) {

	options(shiny.maxRequestSize=1000*1024^2) 

	listObject <- reactiveValues(
		raw=NULL,
		flow.frame = NULL,
		params = NULL,
		id.list=NULL,
		contributions= NULL
	)

	observeEvent(input$upload_file,{
		progress <- Progress$new()
		progress$set(message="Read FCS files...", value=1)
		listObject$raw <- read.FCS(input$upload_file$datapath, emptyValue=FALSE)
		listObject$params <- pData(parameters(listObject$raw))[,1]
		params.names <- pData(parameters(listObject$raw))[,2]
		params.names[which(is.na(params.names))] <- colnames(listObject$raw)[c(which(is.na(params.names)))]
		names(listObject$params) <- params.names
		progress$close()

		output$randDownsampling <-renderUI({
			selectInput("randomDownSampling","Random Downsampling",choices=c(1000,10000,20000,50000,100000))
		})

		output$transformMeth <- renderUI({
			selectInput("transMeth","Select Transform Methode", choices=c("none","arcsinh","logicl"))
		})

		output$markerTransform <- renderUI({
			selectInput("markerTrans","Select Markers to Transform", choices=c("",listObject$params),multiple=TRUE)
		})

	})

	observe({
		if(is.null(input$transMeth))return(NULL)
		if(input$transMeth=="none"){
			disable("markerTrans")
			disable("trans")
		} else if(input$transMeth=="logicl"){
			disable("markerTrans")
			enable("trans")
		} else {
			enable("markerTrans")
			enable("trans")
		}
	})

	observeEvent(input$trans,{
		if(is.null(listObject$flow.frame) || input$transMeth=="none") return(NULL)
		if(is.null(input$markerTrans) && input$transMeth!="logicl") return(NULL)
		progress <- Progress$new()
		progress$set(message="Transform in progress...", value=1)
		if(input$transMeth=="arcsinh"){
			asinhTrans <- arcsinhTransform(transformationId="ln-transformation", a=1, b=1, c=1)
			translist <- transformList(input$markerTrans, asinhTrans) 
			listObject$flow.frame <- transform(listObject$flow.frame, translist)
		} else if(input$transMeth=="logicl") {
			# lgcl <- logicleTransform( w = 0.5, t= 262144, m =4.5)
			# lgcl <- estimateLogicle(listObject$flow.frame, channels =input$markerTrans)
			# trans <- transformList(input$markerTrans, lgcl)
			# listObject$flow.frame <- transform(listObject$flow.frame, lgcl)
			listObject$flow.frame <- logiclTransformCiphe(listObject$flow.frame)
		}
		progress$close()
	})

	observe({
		if(is.null(listObject$raw) || is.null(input$randomDownSampling)) return(NULL)
		progress <- Progress$new()
		progress$set(message="Random Downsampling in progress..", value=1)
		val <- as.integer(input$randomDownSampling)
		val2 <- as.vector(dim(listObject$raw))[1]
		if(val>val2){val <- (val2-1)}
		listObject$flow.frame <- listObject$raw[sample(c(1:dim(listObject$raw)[1]), val),]
		progress$close()
	})

	output$selectID <- renderUI({
		if(is.null(listObject$params)) return(NULL)
		selectInput("id","Select Annotation",choices=c("",listObject$params),multiple=FALSE)
	})

	output$selectIdView <- renderText({
		if(is.null(input$id) || input$id=="") return(NULL)
		listObject$id.list <- unique(listObject$flow.frame@exprs[,input$id])
		return(unique(listObject$flow.frame@exprs[,input$id]))
	})

	output$selectIdValue <- renderUI({
		if(is.null(input$id) || input$id=="") return(NULL)
		selectInput("idValue","Select Id annotation",choices=listObject$id.list)
	})

	output$selectMarkers <- renderUI({
		if(is.null(listObject$params)) return(NULL)
		selectInput("markers","Select Markers",choices=listObject$params,multiple=TRUE)
	})

	observeEvent(input$addAllMarkers,{
		if(is.null(listObject$params)) return(NULL)
		updateSelectInput(session, "markers","Select Markers",
			choices=listObject$params,
			selected=listObject$params
		)
	})

	observeEvent(input$addSameMarkers,{
		if(is.null(listObject$markerTrans)) return(NULL)
		updateSelectInput(session, "markers","Select Markers",
			choices=listObject$params,
			selected=input$markerTrans
		)
	})

	output$selectViewX <- renderUI({
		if(is.null(listObject$params)) return(NULL)
		selectInput("viewX", "X label to view", choices=listObject$params, multiple=FALSE)
	})

	output$selectViewY <- renderUI({
		if(is.null(listObject$params)) return(NULL)
		selectInput("viewY", "Y label to view", choices=listObject$params, multiple=FALSE)
	})

	output$plotPreView <- renderPlot({
		if(is.null(input$viewX) || is.null(input$viewY) || is.null(input$idValue)) return(NULL)
		progress <- Progress$new()
		progress$set(message="Plot in progress...", value=1)
		plot(listObject$flow.frame@exprs[,c(input$viewX, input$viewY)], pch=".",
			col=ifelse(listObject$flow.frame@exprs[,input$id]==input$idValue,"red","blue")
		)
		progress$close()
	},width=400, height=400)

	observeEvent(input$runHypergate,{
		if(is.null(input$idValue) || is.null(input$markers)) return(NULL)

		progress <- Progress$new()
		progress$set(message="Hypergate in progress...", value=1)

		gate_vector <-  rep(0, dim(listObject$flow.frame)[1])
		gate_vector[which(listObject$flow.frame@exprs[,input$id]==input$idValue)] <- 1

		hg_output <- hypergate(
			xp = listObject$flow.frame@exprs[, input$markers],
			gate_vector = gate_vector, 
    	level = 1, verbose = FALSE
    )
    listObject$hg_output

    gating_predicted = subset_matrix_hg(hg_output, listObject$flow.frame@exprs[, input$markers])
  
  	output$tableOutput <- renderTable({
  		table(ifelse(gating_predicted, "Gated-in", "Gated-out"), 
      	ifelse(gate_vector == 1, "Events of interest", "Others")
    	)
    },colnames=FALSE)

  	output$gatingStrat <- renderUI({
  		l <- length(hg_output$active_channels)
  		
  		channels<-hg_output$active_channels
    	channels<-sub("_max","",channels)
    	channels<-sub("_min","",channels)
    	ranges.global<-apply(listObject$flow.frame@exprs[,channels,drop=F],2,range)
    	rownames(ranges.global)=c("min","max")
    	
    	active_events<-rep(T,nrow(listObject$flow.frame@exprs))
    	
    	highlight <- "red"
    	cols<-rep("black",nrow(listObject$flow.frame@exprs))
    	cols[gate_vector==1]=highlight

    	parameters<-hg_output$pars.history
	    active_parameters<-hg_output$active_channels##apply(parameters,2,function(x){x[length(x)]!=x[1]})
	    parameters<-parameters[,active_parameters,drop=FALSE]
	    parameters_order<-apply(parameters,2,function(x)min(which(x!=x[1])))
	    parameters<-parameters[,order(parameters_order,decreasing=FALSE),drop=FALSE]
	    parameters<-setNames(parameters[nrow(parameters),,drop=TRUE],colnames(parameters))

    	channels <-sub("_max","",names(parameters))
    	channels <- sub("_min","",channels)

    	plotOutputObject <- list()

    	n<-length(parameters)
			iter<-0

    	direction <- rep(2,length(parameters))
    	direction[grep("_max",names(parameters))]=1

  		plotOutputObject <- lapply(c(seq(1,n,by=2)), function(i) {

  			tmp <- active_events
  	
  			if((i+1)<=n){
  				iter<<-iter+1
  				chan1<-channels[i]
          chan2<-channels[i+1]

          plotname <- paste("plotGate",iter,sep="")
          plot_output_object <- plotOutput(plotname)
  				plot_output_object <- renderPlot({

  					plot(
  						listObject$flow.frame@exprs[which(tmp==TRUE),chan1],
  						listObject$flow.frame@exprs[which(tmp==TRUE),chan2],
  						xlab=names(listObject$params)[grep(chan1,listObject$params)],
              ylab=names(listObject$params)[grep(chan2,listObject$params)],
              xlim=ranges.global[,chan1],
              ylim=ranges.global[,chan2],
              bty="l",
              pch=20,
              cex=0.1,
              col=cols[tmp],
              main=""
  					)
  					segments(
                x0=parameters[i],
                y0=parameters[i+1],
                x1=ranges.global[direction[i],chan1],
                col="red"
            )
            segments(
                x0=parameters[i],
                y0=parameters[i+1],
                y1=ranges.global[direction[i+1],chan2],
                col="red"
            )
          }, outputArgs = list(width="250px", height ="250px"))

          ## Updating active_events
          if(direction[i]==2){
	          test1<-listObject$flow.frame@exprs[,chan1]>=parameters[i] ##If _min, events above parameter are selected
	        } else {
	          test1<-listObject$flow.frame@exprs[,chan1]<=parameters[i] ##Else events above parameter below
	        }
	        if(direction[i+1]==2){
	          test2<-listObject$flow.frame@exprs[,chan2]>=parameters[i+1]
	        } else {
	          test2<-listObject$flow.frame@exprs[,chan2]<=parameters[i+1]
	        }
	        active_events<<-active_events&test1&test2

	        return(plot_output_object)
  			}

			if(length(parameters)%%2==1){ 

		    chan1<-channels[i]
	      iter<<-iter+1

	      if(direction[iter]==2){
		      test1=listObject$flow.frame@exprs[,chan1]>=parameters[i] ##If _min, events above parameter are selected
		    } else {
		      test1=listObject$flow.frame@exprs[,chan1]<=parameters[i] ##Else events above parameter below
		    }
		    active_events<<-active_events&test1

		     	plotname <- paste("plotGate",iter,sep="")
          plot_output_object <- plotOutput(plotname)
  				plot_output_object <- renderPlot({
			    plot(listObject$flow.frame@exprs[which(active_events==TRUE),chan1],main="",
			      ylab=names(listObject$params)[grep(chan1,listObject$params)],
			      xlab="Events index",
			      ylim=ranges.global[,chan1],
			      bty="l",
			      pch=16,
			      cex=0.1,
			      col=cols[active_events]
			    )
			  	}, outputArgs = list(width="250px", height ="250px"))
		    return(plot_output_object)
		  }

		 })

  		do.call(tagList, plotOutputObject)
  		return(plotOutputObject)
  	})

		contributions <- channels_contributions(hg_output, 
			xp = listObject$flow.frame@exprs[,input$markers], 
			gate_vector = gate_vector, 
    	level = 1, beta = 1
    )

		listObject$contributions <- contributions

       output$barPlot <- renderPlot({
    	channels <- names(contributions) 
    	channels <- sub("_max","",channels)
    	channels <- sub("_min","",channels)
			barplot(contributions,horiz = TRUE,axisnames=F, las=3, space=0,
			ylab = "",
	    xlab = "F1-score deterioration when the parameter is ignored")
	    for(i in 1:length(contributions)){
	    	text(max(contributions)*0.8,(i-(0.1*(1/i))),names(listObject$params)[which(listObject$params==channels[i])],cex=1.25)
	    }
		},height=100*length(hg_output$active_channels))
		progress$close()

		if(!is.null(dim(hg_output$pars.history[,hg_output$active_channels]))){
			output$optiButton <- renderUI({
				actionButton("opti","Reoptimize")
			})
		}
	})

	observeEvent(input$opti,{
		hg_output <- listObject$hg_output
		contributions <- listObject$contributions
		gate_vector <- rep(0, dim(listObject$flow.frame)[1])
		gate_vector[which(listObject$flow.frame@exprs[,input$id]==input$idValue)] <- 1

		 hg_output <- reoptimize_strategy(gate = hg_output, 
          channels_subset = names(contributions)[c(1,2)], 
          xp = listObject$flow.frame@exprs[,input$markers],
          gate_vector = gate_vector, level = 1)

		output$gatingStrat <- renderUI({
  		l <- length(hg_output$active_channels)
  		iter<-0
  		channels<-hg_output$active_channels
    	channels<-sub("_max","",channels)
    	channels<-sub("_min","",channels)
    	ranges.global<-apply(listObject$flow.frame@exprs[,channels,drop=F],2,range)
    	rownames(ranges.global)=c("min","max")
    	active_events<-rep(T,nrow(listObject$flow.frame@exprs))
    	highlight <- "red"
    	cols<-rep("black",nrow(listObject$flow.frame@exprs))
    	cols[gate_vector==1]=highlight

  		plotOutputObject <- lapply(c(1:l), function(i){

  			if((i+1)<=l){
  				chan1<-channels[i]
          chan2<-channels[i+1]
          iter<-iter+1

  				renderPlot({
  					plot(
  						listObject$flow.frame@exprs[active_events,chan1],
  						listObject$flow.frame@exprs[active_events,chan2],
  						xlab=names(listObject$params)[grep(chan1,listObject$params)],
              ylab=names(listObject$params)[grep(chan2,listObject$params)],
              xlim=ranges.global[,chan1],
              ylim=ranges.global[,chan2],
              bty="l",
              pch=20,
              cex=0.1,
              col=cols[active_events]
  					)
          }, outputArgs = list(width="250px", height ="250px"))
  			}
  		})
  		do.call(tagList, plotOutputObject)
  		return(plotOutputObject)
  	})
	})

}