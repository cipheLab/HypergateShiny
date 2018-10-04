library("shinyjs")

ui <- fluidPage(theme = "hypergate.css",
	shinyjs::useShinyjs(),
	titlePanel("Hypergate"),
	sidebarLayout(
		sidebarPanel(width = 3,
			fileInput("upload_file","Upload File (FCS)",
        multiple = FALSE,
        accept = c("text/fcs",".fcs"),
        placeholder = "no file selected"
      ),
      uiOutput("randDownsampling"),
      uiOutput("transformMeth"),
      uiOutput("markerTransform"),
      actionButton("trans","Transform"),
      uiOutput("selectID"),
      textOutput("selectIdView"),
      uiOutput("selectIdValue"),
      uiOutput("selectMarkers"),
      actionButton("addAllMarkers","Add All"),
      actionButton("addSameMarkers","Add Same Markers"),
      tags$br(),tags$br(),
      uiOutput("selectViewX"),
      uiOutput("selectViewY")
		),
		mainPanel(
			fluidRow(actionButton("runHypergate","Run Hypergate")),
			fluidRow(
				column(5,
					fluidRow(plotOutput("plotPreView",height="400px")),
					fluidRow(tableOutput("tableOutput")),
					fluidRow(uiOutput("gatingStrat")),
					fluidRow(uiOutput("optiButton"))
				),
				column(7,
					fluidRow(plotOutput("barPlot",height="100%"))
				)
			)
		)
	)
)