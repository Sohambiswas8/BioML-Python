install.packages('EBImage')
library(BiocManager)
BiocManager::install("EBImage")
library(EBImage)

# Contrast function
contrast <- function(image,f){
  imageData(image)[,,1] = ((imageData(image)[,,1] - 0.5) * (1 + f)) + 0.5 
  imageData(image)[,,2] = ((imageData(image)[,,2] - 0.5) * (1 + f)) + 0.5 
  imageData(image)[,,3] = ((imageData(image)[,,3] - 0.5) * (1 + f)) + 0.5 
  return(image)
}
# brightness
bright = function(image,b){
  imageData(image)[,,1] =  ifelse( (imageData(image)[,,1]+b) < 0,0,ifelse( (imageData(image)[,,1]+b) > 1 ,1, (imageData(image)[,,1]+b)) )
  imageData(image)[,,2] =  ifelse( (imageData(image)[,,2]+b) < 0,0,ifelse( (imageData(image)[,,2]+b) > 1 ,1, (imageData(image)[,,2]+b)) )
  imageData(image)[,,3] =  ifelse( (imageData(image)[,,3]+b) < 0,0,ifelse( (imageData(image)[,,3]+b) > 1 ,1, (imageData(image)[,,3]+b)) )
  return(image)
}
# gray scale 
gray = function(image){
  newMat = ((imageData(image)[,,1]+imageData(image)[,,2]+imageData(image)[,,3])/3)
  newImg = EBImage::Image(newMat, c(nrow(newMat),ncol(newMat),1), colormode = "Gray")
  return(newImg)
}
# red Channel : making green and blue 'zero'
red = function(image){
  imageData(image)[,,2] =  0
  imageData(image)[,,3] =  0
  return(image)
}
# green Channel : making red and blue 'zero'
green = function(image){
  imageData(image)[,,1] =  0
  imageData(image)[,,3] =  0
  return(image)
}
# blue Channel : making red and green 'zero'
blue = function(image){
  imageData(image)[,,1] =  0
  imageData(image)[,,2] =  0
  return(image)
}


ui <- fluidPage(
  titlePanel("Manipulate images using pixel values in R"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("contrast_factor", 
                  "Contrast Adjustment:", 
                  min = -1, 
                  max = 1, 
                  value = 0, 
                  step = 0.1),
      sliderInput("brightness_factor", 
                  "Brightness Adjustment:", 
                  min = -1, 
                  max = 1, 
                  value = 0, 
                  step = 0.1),
      radioButtons("channel", "Choose channel:", 
                   choices = c("All" = "all",
                               "Red" = "red", 
                               "Green" = "green", 
                               "Blue" = "blue",
                               "Make Gray Scale" = "gray")
      )
    ),
    mainPanel(
      plotOutput("nebula_plot")
    )
  )
)

server <- function(input, output) {
  
  # Read image
  nebula <- EBImage::readImage('https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/Crab_Nebula.jpg/1920px-Crab_Nebula.jpg')
  
  # Render plot output
  output$nebula_plot <- renderPlot({
    image <- nebula
    image <- contrast(image, input$contrast_factor)
    image <- bright(image, input$brightness_factor)
    if(input$channel == "red") {
      image<-red(image)
    }else if(input$channel == "green") {
      image<-green(image)
    }else if(input$channel == "blue") {
      image<-blue(image)
    }else if(input$channel == "gray") {
      image<-gray(image)
    }else{
      image<-image
    }
    
    plot(image)
  })
}

shiny::shinyApp(ui = ui, server = server)

