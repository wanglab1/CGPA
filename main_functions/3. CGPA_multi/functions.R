###################### Useful functions ######################################
####### Func1: navdash with search box ########

navbarPageWithInputs <- function(..., inputs) {
  navbar <- navbarPage(...)
  form <- tags$form(class = "navbar-form", inputs)
  navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
    navbar[[3]][[1]]$children[[1]], form)
  navbar
}




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
