###################### Useful functions ######################################
####### Func1: navdash with search box ########
############ Function to include searchbox on top ######################
# navbarPageWithInputs <- function(..., inputs1,inputs2) {
#   navbar <- navbarPage(...)
#   form1 <- tags$form(class = "navbar-form", inputs1)
#   navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
#     navbar[[3]][[1]]$children[[1]], form1)
#   form2 <- tags$form(class = "navbar-form", inputs2)
#   navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
#     navbar[[3]][[1]]$children[[1]], form2)
#   navbar
# }

navbarPageWithInputs <- function(..., inputs) {
  navbar <- navbarPage(...)
  form <- tags$form(class = "navbar-form", inputs)
  navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
    navbar[[3]][[1]]$children[[1]], form)
  navbar
}

###### func2: dark background image ############
library(ggplot2)
library(gridExtra)

# # For forest plot
# theme_black = function(base_size = 12, base_family = "") {
#   
#   theme_grey(base_size = base_size, base_family = base_family) %+replace%
#     
#     theme(
#       # Specify axis options
#       axis.line = element_line(colour = "white",linetype = "solid"),  
#       axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
#       axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9,face="bold"),  
#       axis.ticks = element_line(color = "white", size  =  0.2),  
#       axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
#       axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
#       axis.ticks.length = unit(0.3, "lines"),   
#       # Specify legend options
#       legend.background = element_rect(color = NA, fill = "black"),  
#       legend.key = element_rect(color = "white",  fill = "black"),  
#       legend.key.size = unit(1.2, "lines"),  
#       legend.key.height = NULL,  
#       legend.key.width = NULL,      
#       legend.text = element_text(size = base_size*0.8, color = "white"),  
#       legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
#       legend.position = "none",  
#       legend.text.align = NULL,  
#       legend.title.align = NULL,  
#       legend.direction = "vertical",  
#       legend.box = NULL, 
#       # Specify panel options
#        panel.background = element_rect(fill = "black", color  =  NA),  
#       # panel.border = element_rect(fill = NA, color = "white"),  
#       # panel.grid.major = element_line(color = "grey35"),  
#       # panel.grid.minor = element_line(color = "grey20"),  
#       panel.border = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.margin = unit(0.5, "lines"),   
#   
#       # Specify facetting options
#       strip.background = element_rect(fill = "grey30", color = "grey10"),  
#       strip.text.x = element_text(size = base_size*0.8, color = "white"),  
#       strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
#       # Specify plot options
#       plot.background = element_rect(color = "black", fill = "black"),  
#       plot.title = element_text(size = base_size*1.2, color = "white"),  
#       plot.margin = unit(rep(1, 4), "lines")
#       
#     )
#   
# }
# 
# # for heatplot
# theme_black2 = function(base_size = 12, base_family = "") {
#   
#   theme_grey(base_size = base_size, base_family = base_family) %+replace%
# 
#     theme(
#       # Specify axis options
#       axis.line = element_line(colour = "white",linetype = "solid"),  
#       axis.title.y=element_blank(),
#       axis.title.x = element_text(size = base_size, color = "black", margin = margin(0, 10, 0, 0)),  
#       axis.text.x = element_text(size = base_size*0.8, color = "black", lineheight = 0.9),  
#       axis.ticks.x=element_blank(),
#       axis.text.y = element_text(face="bold", color="white"),
# 
#     
#       # Specify legend options
#       legend.background = element_rect(color = NA, fill = "black"),  
#       legend.key = element_rect(color = "white",  fill = "black"),  
#       legend.key.size = unit(1.2, "lines"),  
#       legend.key.height = NULL,  
#       legend.key.width = NULL,      
#       legend.text = element_text(size = base_size*0.8, color = "white"),  
#       legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
#     #  legend.position = "none",  
#       legend.text.align = NULL,  
#       legend.title.align = NULL,  
#       legend.direction = "vertical",  
#       legend.box = NULL, 
#       # Specify panel options
#       panel.background = element_rect(fill = "black", color  =  NA),  
#       # panel.border = element_rect(fill = NA, color = "white"),  
#       # panel.grid.major = element_line(color = "grey35"),  
#       # panel.grid.minor = element_line(color = "grey20"),  
#       panel.border = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.margin = unit(0.5, "lines"),   
#       
#       # Specify facetting options
#       strip.background = element_rect(fill = "grey30", color = "grey10"),  
#       strip.text.x = element_text(size = base_size*0.8, color = "white"),  
#       strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
#       # Specify plot options
#       plot.background = element_rect(color = "black", fill = "black"),  
#       plot.title = element_text(size = base_size*1.2, color = "white"),  
#       plot.margin = unit(rep(1, 4), "lines")
#       
#     )
#   
# }
# 
# # for bodymap 
# theme_black3 = function(base_size = 12, base_family = "") {
#   
#   theme_grey(base_size = base_size, base_family = base_family) %+replace%
#     
#     theme(
#       # Specify axis options
#       axis.line = element_blank(),  
#       axis.title.y=element_blank(),
#       axis.title.x = element_blank(),
#       axis.text.x = element_blank(),  
#       axis.ticks.x=element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y=element_blank(),
#       
#       
#       # Specify legend options
#       
#       legend.background = element_rect(color = NA, fill = "#696969"),  
#       legend.key = element_rect(color = "white",  fill = "#696969"),  
#       legend.key.size = unit(1.2, "lines"),  
#       legend.key.height = NULL,  
#       legend.key.width = NULL,      
#       legend.text = element_text(size = base_size*0.8, color = "white"),  
#       legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
#        legend.position = "none",  
#       legend.text.align = NULL,  
#       legend.title.align = NULL,  
#       legend.direction = "vertical",  
#       legend.box = NULL, 
#       # Specify panel options
#       panel.background = element_rect(fill = "transparent",colour = NA),
#       # panel.border = element_rect(fill = NA, color = "white"),  
#       # panel.grid.major = element_line(color = "grey35"),  
#       # panel.grid.minor = element_line(color = "grey20"),  
#       panel.border = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
# 
#       # Specify plot options
#       plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
# 
#       
#     )
#   
# }
#  
# 
# # for cir bar 
# theme_black4 = function(base_size = 12, base_family = "") {
#   
#   theme_grey(base_size = base_size, base_family = base_family) %+replace%
#     
#     theme(
#       # Specify axis options
#       axis.line = element_blank(),  
#       axis.title.y=element_blank(),
#       axis.title.x = element_blank(),
#       axis.text.x = element_blank(),  
#       axis.ticks.x=element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y=element_blank(),
#       
#       
#       # Specify legend options
#       
#       legend.background = element_rect(fill = "transparent"),
#       legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#       
#       legend.key = element_rect(color = "white",  fill = "#696969"),  
#       legend.key.size = unit(1.2, "lines"),  
#       legend.key.height = NULL,  
#       legend.key.width = NULL,      
#       legend.text = element_text(size = base_size*0.8, color = "white"),  
#       legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
#       legend.position = "none",  
#       legend.text.align = NULL,  
#       legend.title.align = NULL,  
#       legend.direction = "vertical",  
#       legend.box = NULL, 
#       # Specify panel options
#       panel.background = element_rect(fill = "transparent",colour = NA),
#       # panel.border = element_rect(fill = NA, color = "white"),  
#       # panel.grid.major = element_line(color = "grey35"),  
#       # panel.grid.minor = element_line(color = "grey20"),  
#       panel.border = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#     #  plot.margin = unit(rep(-1,4), "cm"),
#       # Specify plot options
#       plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      
#       
#     )
#   
# }

#===============================================================================
######################## common functions for UI ##########################
#===============================================================================
########### functions for selecting genes ###############
sel_input_ui = function(input_id,input_name,choices = NULL,multiple=FALSE,selected = NULL){
  selectizeInput(input_id,input_name,choices=choices,multiple=multiple,selected=selected)
}

sel_input_server = function(input,output,session,input_id,input_name,input_res,sel_input){
  updateSelectizeInput(session,input_id,input_name,choices=input_res,selected=sel_input,server = T)
  
}
