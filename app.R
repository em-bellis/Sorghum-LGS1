#######This is a demo R shiny app based on a species distribution model for the parasitic weed Striga hermonthica.
#######The two dependent files (preds.all.tif and curated_occ_ESB_7.26.18.csv) can be downloaded from my github page: 
#######Last update: 4.2.19 by Emily Bellis

#######to deploy in test mode from R:
#library(shiny)
#runApp("/Users/emilywork/Desktop/all_occ_model")
#######after each change can save app.R and just reload browser window to view changes
#######for more info, see https://shiny.rstudio.com/articles/build.html

library(shiny)
library(shinythemes)
library(raster)
library(maps)

all <- raster('/Users/emilywork/Desktop/ENM/distModel7.30/preds.all.tif') ##raster layer output from species distribution model
sorg <- raster('/Users/emilywork/Desktop/ENM/distModel7.30/preds.sorg.tif', header=T)
shgeo <- read.csv('/Users/emilywork/Desktop/ENM/Final_occurrence_data/curated_occ_ESB_7.26.18.csv', header=T,stringsAsFactors=F) ##1369 S. hermonthica occurrences
data(worldMapEnv)

# Define UI/layout ----
ui <- fluidPage(theme=shinytheme("spacelab"),
  # App title ----
  titlePanel(HTML("<em>Striga hermonthica</em> distribution model")),
 
  # Sidebar panel for inputs and outputs----
  sidebarPanel(
     # Three types of inputs ----
      numericInput(inputId = "lat",
                  label = "Latitude (decimal degrees)",
                  value = -0.10),
      numericInput(inputId = "lon",
                  label = "Longitude (decimal degrees)",
                  value = 34.76),
      checkboxInput(inputId = "occs", "Plot known occurrence records?", value=FALSE
      		),
      checkboxInput(inputId = "sorgmod", "Use sorghum-only model?", value=FALSE
      		),
      		
      # Output panel
      h4("Habitat suitability score at location of interest is", textOutput("hs"))
  ),

  # Main panel for displaying outputs ----
  mainPanel(
  	plotOutput("baseplot"),
    p(("Maxent species distribution model based on") ,em("S. hermonthica "),("occurrences and 8 environmental variables (annual rainfall, soil nitrogen, soil phosphorus, soil clay fraction, mean temperature of the wettest quarter, isothermality, topographic wetness index, and PET). The default model is based on 1050 parasite occurrence records. The sorghum-only model excludes parasites on unknown or non-sorghum hosts, leaving a subset of 262 parasite occurrence records.")),
  	p("Additional information about the models can be found in Bellis et al. (in prep)."),
  	p("Plot may take a minute or so to load!")
  )
)

# Define server logic to plot things with reactivity
server <- function(input, output) {
	#put all the reactive expressions up here
		getbasemap <- reactive({ #this is in a reactive expression because it is shared by both the plot and function that extracts the value of the model
			if(input$sorgmod==TRUE){  #this part responds to the check box if user wants the sorghum-only model
				sorg
			} else {all}
		})
		Lat <- reactive({input$lat})
		Lon <- reactive({input$lon})
		Lat_d <- debounce(Lat, 2000) #slow this down to 3 seconds allow user to type
		Lon_d <- debounce(Lon, 2000) #slow this down to allow user to type
		
	#the first render thing I want to do is my main plot
	output$baseplot <- renderPlot({
		map(database="world",lwd=0.5, col="grey50")
		basemap <- getbasemap()
		plot(basemap, add=T, legend=F)
		map(database="world", add=T, lwd=0.5, col="grey50")
		
		if(input$occs==TRUE){  #this part responds to the check box if user wants to plot all the S. hermonthica occurrences
			points(shgeo$lon, shgeo$lat, col='black',cex=0.1, pch=20)
		}
		
		points(x=Lon_d(), y=Lat_d(), col='black', pch='X') #mark location in side panel with an X
		scalebar(5000, xy = c(100,-57), type = "bar", divs = 2, lonlat = TRUE)
		plot(basemap, legend.only=T, add=T, smallplot=c(0.2,0.25, 0.3,0.45))
		text(-120,-7, "Habitat suitability")
		text(100,-62, "km")
	})
	
	#the second render thing is extracting the value of the model at any gps location
	output$hs <- renderText({
		basemap <- getbasemap()
		extract(basemap,cbind(Lon_d(),Lat_d()))
	})
}

shinyApp(ui, server)