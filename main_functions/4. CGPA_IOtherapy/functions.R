###################### Useful functions ######################################
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

sel_input_server_fixed = function(input,output,session,input_id,input_name,input_res,sel_input){
  updateSelectizeInput(session,input_id,input_name,choices=input_res,selected=sel_input,server = T,options = list(maxItems  = 5))
  
}
# updateSelectize for renderUI
update_ui_select = function(session,input_info,geneID){
  updateSelectizeInput(session,input_info,choices = geneID,server = T)
  
}

# Define the mapping function to rename study names
map_study_names <- function(study_name) {
  study_name <- gsub(" \\(", "_", study_name)  # Replace " (" with "_"
  study_name <- gsub("\\)", "", study_name)    # Remove ")"
  study_name
}